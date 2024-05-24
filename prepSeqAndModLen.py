import json
import argparse
import sys
import os

# Import SQLAlchemy
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

# Import our schemas
sys.path.append(os.path.join(os.path.dirname(__file__), "../Schemata/ORMs/python"))
import dfamorm as dfam
import assemblydborm as adbORM

# Import Dfam Libraries
sys.path.append(os.path.join(os.path.dirname(__file__), "../Lib"))
import DfamConfig as dc
import DfamVersion as dfVersion
import DfamDBView as dfv


def _usage():
    """Print out docstring as program usage"""
    # Call to help/pydoc with scriptname ( sans path and file extension )
    help(os.path.splitext(os.path.basename(__file__))[0])
    sys.exit(0)


def main(*args):
    # Options processing
    #
    #   There are two ways to document usage/command line
    #   arguments using this boilerplate.  The intended way
    #   is to document using docstrings at the top of the
    #   script.  This way the pydoc docs match the output
    #   produced by '-h' or '--help' using the argparse
    #   custom action class ( _CustomUsageAction ) defined
    #   below.  If you want the paired-down argparse default
    #   instead simply remove the "add_help=False" argument
    #   to the argparse constructor below and comment out
    #   the add_argment('-h', ...) line below.
    #
    class _CustomUsageAction(argparse.Action):
        def __init__(
            self, option_strings, dest, default=False, required=False, help=None
        ):
            super(_CustomUsageAction, self).__init__(
                option_strings=option_strings,
                dest=dest,
                nargs=0,
                const=True,
                default=default,
                required=required,
                help=help,
            )

        def __call__(self, parser, args, values, option_string=None):
            _usage()

    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("-h", "--help", action=_CustomUsageAction)
    parser.add_argument("-c", "--dfam_config", dest="dfam_config")
    parser.add_argument("-a", "--assembly", dest="assembly", required=True)
    parser.add_argument("-r", "--dfam_release", dest="dfam_release", required=True)
    parser.add_argument("-d", "--debug", dest="debug", action="store_true")
    parser.add_argument("-q", "--quiet", dest="quiet", action="store_true")
    parser.add_argument("-v", "--version", dest="get_version", action="store_true")
    #
    # Examples:
    #   e.g. -f 3
    #     parser.add_argument('-f','--foo', type=int, default=42, help='FOO!')
    #   e.g. -f   : If set store True in args.foo
    #     parser.add_argument('-f','--foo', action='store_true')
    #   e.g. Set the -foo parameter as required
    #     parser.add_argument('-f','--foo', type=int, required=True )
    #   e.g. Set the args.foobar using the parameter passed with -foo
    #     parser.add_argument('-f','--foo', dest='foobar', type=int )
    #
    args = parser.parse_args()

    #
    # Require Dfam Configuration
    #   Search order: --dfam_config <path>, environment DFAM_CONF,
    #                 and finally "../Conf/dfam.conf"
    #
    conf = dc.DfamConfig(args.dfam_config)
    df_ver = dfVersion.DfamVersion()

    if not os.path.exists(".my.cnf"):
        with open(".my.cnf", "w") as cnf:
            cnf.write("[mysql]\n")
            cnf.write("user=" + conf.getKeyValue("schema.Dfam.user") + "\n")
            cnf.write("password=" + conf.getKeyValue("schema.Dfam.password") + "\n")
        os.chmod(".my.cnf", 0o600)

    if args.get_version:
        print(df_ver.version_string)
        exit(0)

    schema_suffix = "_df" + args.dfam_release.replace(".", "")

    #
    # Announce ourselves unless otherwise requested
    #
    if not args.quiet:
        print("#\n# buildFullRegion.py " + df_ver.version_string + "\n#")
        print("Parameters:")
        print(" - assembly: " + args.assembly)
        print(" - dfam release: " + args.dfam_release)
        print(" - sql schema target: " + args.assembly + schema_suffix)

    #
    # Database Connection ( required )
    #
    df_engine = create_engine(conf.getDBConnStrWPassFallback("Dfam"))
    df_sfactory = sessionmaker(df_engine)
    df_session = df_sfactory()

    # Obtain assembly and taxonomy starting point
    assembly = (
        df_session.query(dfam.Assembly)
        .filter(dfam.Assembly.name == args.assembly)
        .one()
    )
    dateformat = "%y-%m-%d %H:%M:%S"
    allfams = [fam[1] for fam in dfv.getFamiliesInAssembly(df_session, args.assembly)]
    print("Families identified in " + args.assembly + " : " + str(len(allfams)))

    mod_len_info = {
        "assembly": assembly.schema_name,
        "version": assembly.version,
        "release_date": assembly.release_date.strftime(dateformat),
    }
    mod_len_json_file = f"{args.assembly}-model_lengths.json"

    seq_info = {
        "assembly": assembly.schema_name,
        "version": assembly.version,
        "release_date": assembly.release_date.strftime(dateformat),
    }
    seq_json_file = f"{args.assembly}-sequences.json"

    print(f"Extracting Model Length Data - {args.assembly}")

    mod_len_data = {
        fam[0]: {"length": fam[1]}
        for fam in df_session.query(dfam.Family.accession, dfam.Family.length)
        .filter(dfam.Family.accession.in_(allfams))
        .all()
    }
    mod_len_info["data"] = mod_len_data

    print(f"Extracting Sequence Data - {args.assembly}")
    adb_engine = create_engine(conf.getADBConnStrWPassFallback(assembly.schema_name))
    adb_sfactory = sessionmaker(adb_engine)
    adb_session = adb_sfactory()
    aseqs = adb_session.query(adbORM.Sequence).all()

    seq_data = {
        seq.accession: {
            "id": seq.id,
            "description": seq.description,
            "length": seq.length,
            "updated": seq.updated.strftime(dateformat),
            "created": seq.created.strftime(dateformat),
            "is_genomic": seq.is_genomic,
        }
        for seq in aseqs
    }
    seq_info["data"] = seq_data

    print(f"Saving Sequence JSON - {seq_json_file}")
    with open(seq_json_file, "w") as fp:
        json.dump(seq_info, fp)

    print(f"Saving Model Length JSON - {mod_len_json_file}")
    with open(mod_len_json_file, "w") as fp:
        json.dump(mod_len_info, fp)


#
# Wrap script functionality in main() to avoid automatic execution
# when imported ( e.g. when help is called on file )
#
if __name__ == "__main__":
    main(*sys.argv)
