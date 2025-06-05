"""
Command Line Interface
"""
import sys
from argparse import ArgumentParser
import kpfcc

def main():
    # This parser defines the top-level command "astroq"
    psr = ArgumentParser(
        description='AstroQ: Optimized observation scheduling',
        prog='astroq'
    )

    psr.add_argument('-V', '--version',
                     action='version',
                     version="%(prog)s {}".format(kpfcc.__version__),
                     help="Print version number and exit."
                     )
    # Parsers spawned from this subparser generator will define subcommands of the astroq top-level command
    # Parsers spawned from subpsr (so it defines a subcommand) has psr_parent as its parent (so it accepts any arguments of psr_parent in addition to its own)
    subpsr = psr.add_subparsers(title='subcommands', dest='subcommand', required=True)
    ## Parent parser to define arguments common to all subcommands (so you don't have to redefine them)
    psr_parent = ArgumentParser(add_help=False)

    ## subcommand of astroq: bench -- conduct benchmark tests
    psr_bench = subpsr.add_parser('bench', parents=[psr_parent],
                                  description='Conduct benchmark tests',
                                  prefix_chars='-'
                                  )
    psr_bench.add_argument('-nr', '--number_requests',
                            type=int,
                            required=True,
                            help="Run benchmark with given number of requests (strongly suggest value >= 250.)"
                            )

    psr_bench.add_argument('-ns', '--number_slots',
                            type=int,
                            required=True,
                            help="Run benchmark with given number of slots required to complete the extra single shot requests."
                            )

    psr_bench.add_argument('-sc', '--shortcut',
                            type=int,
                            required=False,
                            help="Run benchmark with a small request set, for testing purposes only."
                            )

    psr_bench.add_argument('-cf', '--config_file',
                              type=str,
                              required=True,
                              help="Relative path of config file."
                            )
    psr_bench.set_defaults(func=kpfcc.driver.bench)

    ## subcommand of astroq: plot -- run the plotting suite
    psr_plot = subpsr.add_parser('plot', parents=[psr_parent],
                                  description='Run the plotting suite',
                                  prefix_chars='-'
                                  )
    psr_plot.add_argument('-cf', '--config_file',
                              type=str,
                              required=True,
                              help="Relative path of config file."
                              )
    psr_plot.set_defaults(func=kpfcc.driver.plot)

    ## subcommand of astroq: schedule -- Schedule observation requests
    psr_schedule = subpsr.add_parser('schedule', parents=[psr_parent],
                                      description="Schedule observation requests",
                                      prefix_chars="-"
                                     )
    psr_schedule.add_argument('-rf', '--request_file',
                              type=str,
                              required=True,
                              help="Relative path of request file."
                              )
    psr_schedule.add_argument('-cf', '--config_file',
                              type=str,
                              required=True,
                              help="Relative path of config file."
                              )
    psr_schedule.set_defaults(func=kpfcc.driver.schedule)


    ## subcommand of astroq: kpfcc -- Do KPFCC stuff
    psr_kpfcc = subpsr.add_parser('kpfcc', parents=[psr_parent],
                                  description='Do KPFCC stuff',
                                  prefix_chars='-'
                                  )
    kpfcc_subpsr = psr_kpfcc.add_subparsers(title='kpfcc subcommands', dest='kpfcc_subcommand')
    psr_kpfcc.set_defaults(func=kpfcc.driver.kpfcc)

    ## subcommand of kpfcc: build -- Build observation requests
    psr_kpfcc_build = kpfcc_subpsr.add_parser('build', #parents=[psr_parent],
                                               description="Build observation requests",
                                               prefix_chars="-"
                                               )
    psr_kpfcc_build.add_argument('-cf', '--config_file',
                              type=str,
                              required=True,
                              help="Relative path of config file."
                                )
    psr_kpfcc_build.set_defaults(func=kpfcc.driver.kpfcc_build)

    ## subcommand of kpfcc: prepare -- Prep for a new semester
    psr_kpfcc_prep = kpfcc_subpsr.add_parser('prep', #parents=[psr_parent],
                                               description="Prepare for a new semester",
                                               prefix_chars="-"
                                               )
    psr_kpfcc_prep.add_argument('-cf', '--config_file',
                              type=str,
                              required=True,
                              help="Relative path of config file."
                                )
    psr_kpfcc_prep.set_defaults(func=kpfcc.driver.kpfcc_prep)

    ## subcommand of kpfcc: data -- pull latest OB database
    psr_kpfcc_data = kpfcc_subpsr.add_parser('data', #parents=[psr_parent],
                                               description="Pull the OB database from Keck",
                                               prefix_chars="-"
                                               )

    psr_kpfcc_data.add_argument('-pf', '--pull_file',
                              type=str,
                              required=True,
                              help="Path to the file that determines how to pull from the database."
                                )

    psr_kpfcc_data.add_argument('-df', '--database_file',
                              type=str,
                              required=True,
                              help="Path to save the good OBs request sheet."
                                )

    psr_kpfcc_data.set_defaults(func=kpfcc.driver.kpfcc_data)
    
    ## subcommand of kpfcc: webapp -- launch web app to view interactive plots
    psr_kpfcc_webapp = kpfcc_subpsr.add_parser('webapp', #parents=[psr_parent],
                                               description="Launch web app to view interactive plots",
                                               prefix_chars="-"
                                               )

    psr_kpfcc_webapp.add_argument('-pp', '--pickle_path',
                              type=str,
                              required=True,
                              help="Path to directory containing star_objects.pkl containing data to plot stellar observing info."
                                )

    psr_kpfcc_webapp.add_argument('-cf', '--config_file',
                              type=str,
                              required=True,
                              help="Path to config file."
                                )

    psr_kpfcc_webapp.set_defaults(func=kpfcc.driver.kpfcc_webapp)

    ## subcommand of astroq: ttp -- run the ttp
    psr_ttp = subpsr.add_parser('ttp', parents=[psr_parent],
                                  description='Run the ttp',
                                  prefix_chars='-'
                                  )
    psr_ttp.add_argument('-cf', '--config_file',
                              type=str,
                              required=True,
                              help="Relative path of config file."
                              )
    psr_ttp.set_defaults(func=kpfcc.driver.ttp)
    
    ## subcommand of astroq: compare -- compare request set and schedule file
    psr_compare = subpsr.add_parser('comp', parents=[psr_parent],
                                    description='Compare request set and schedule for consistency',
                                    prefix_chars='-'
                                    )
    psr_compare.add_argument('-rf', '--request_file',
                              type=str,
                              required=True,
                              help="Relative path of request set file."
                              )
    psr_compare.add_argument('-sf', '--schedule_file',
                              type=str,
                              required=True,
                              help="Relative path of schedule file."
                              )
    psr_compare.set_defaults(func=kpfcc.driver.requests_vs_schedule)

    # If no arguments are provided, print help message and exit
    if len(sys.argv)==1:
        psr.print_help(sys.stderr)
        sys.exit(1)

    args = psr.parse_args()
    args.func(args)

if __name__=="__main__":
    main()
