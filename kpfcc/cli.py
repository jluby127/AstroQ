"""
Command Line Interface
"""

import sys
from argparse import ArgumentParser
import kpfcc

# import pdb; pdb.set_trace()

def main():
    # This parser defines the top-level command "astroq"
    psr = ArgumentParser(
        description='AstroQ: Optimized observation scheduling',
        prog='astroq'
    )

    psr.add_argument('-v', '--version',
                     action='version',
                     version="%(prog)s {}".format(kpfcc.__version__),
                     help="Print version number and exit."
                     )
    # Parsers spawned from this subparser generator will define subcommands of the astroq top-level command
    subpsr = psr.add_subparsers(title='subcommands', dest='subcommand')


    ## Parent parser to define arguments common to all subcommands (so you don't have to redefine them)
    ################################################################
    psr_parent = ArgumentParser(add_help=False)

    ## Example arguments to add for all subcommands
    psr_parent.add_argument('-od', '--outdir',
                            type=str,
                            required=False,
                            default=None,
                            help='Path to directory where output files will be saved'
                            )
    # psr_parent.add_argument('-v', '--verbose',
    #                         action='store_true',
    #                         required=False,
    #                         help='Print out warnings and (possibly) useful information'
    #                         )


    ## subcommand of astroq: bench -- conduct benchmark tests
    ############################################
    # Note: this parser is spawned from subpsr (so it defines a subcommand) and has psr_parent as its parent (so it accepts any arguments of psr_parent in addition to its own)
    psr_bench = subpsr.add_parser('bench', parents=[psr_parent],
                                  description='Conduct benchmark tests',
                                  prefix_chars='-'
                                  )

    psr_bench.add_argument('-bm', '--benchmark',
                            type=str,
                            required=True,
                            help="Conduct the specified benchmark test"
                            )

    psr_bench.set_defaults(func=kpfcc.driver.bench)
    
    ## subcommand of astroq: schedule-request -- Schedule observation requests
    #############################################################
    
    psr_schedule = subpsr.add_parser('schedule', parents=[psr_parent],
                                      description="Schedule observation requests",
                                      prefix_chars="-"
                                                  )
                                   
    psr_schedule.add_argument('-rf', '--request_file',
                              type=int,
                              required=True,
                              help="Relative path of request file."
                                    )
    psr_schedule.add_argument('-cf', '--config_file',
                              type=int,
                              required=True,
                              help="Relative path of config file."
                                    )
                                    
    psr_schedule.set_defaults(func=kpfcc.driver.schedule)
    
    
    ## subcommand of astroq: kpfcc -- Do KPFCC stuff
    ############################################
    # Note: this parser is spawned from subpsr (so it defines a subcommand) and has psr_parent as its parent (so it accepts any arguments of psr_parent in addition to its own)
    psr_kpfcc = subpsr.add_parser('kpfcc', parents=[psr_parent],
                                  description='Do KPFCC stuff',
                                  prefix_chars='-'
                                  )
    ## Add KPF CC arguments here
    # psr_kpfcc.add_argument('-arg', '--argument',
    #                         type=str,
    #                         required=True,
    #                         help="Argument to do KPFCC stuff"
    #                         )


    # Parsers spawned from this subparser will define subcommands of the kpfcc subcommand
    kpfcc_subpsr = psr_kpfcc.add_subparsers(title='kpfcc subcommands', dest='kpfcc_subcommand')


    psr_kpfcc.set_defaults(func=kpfcc.driver.kpfcc)

    ## subcommand of kpfcc: build-request: Build observation requests
    #############################################################

    psr_kpfcc_build = kpfcc_subpsr.add_parser('build-request', #parents=[psr_parent],
                                               description="Build observation requests",
                                               prefix_chars="-"
                                               )

    psr_kpfcc_build.add_argument('-boolarg', '--boolean_argument',
                                action='store_true',
                                required=False,
                                help="Set to True if flag is provided; default to False."
                                )

    psr_kpfcc_build.set_defaults(func=kpfcc.driver.kpfcc_build)


    ## subcommand of kpfcc: plot: Make plots of scheduling results
    #############################################################

    psr_kpfcc_plot = kpfcc_subpsr.add_parser('plot', #parents=[psr_parent],
                                             description="Plot scheduling results",
                                             prefix_chars="-"
                                             )

    psr_kpfcc_plot.add_argument('-so', '--schedule_object',
                                type=str,
                                required=True,
                                help="Whatever object is needed to make plots."
                                )


    psr_kpfcc_plot.set_defaults(func=kpfcc.driver.kpfcc_plot)


    # If no arguments are provided, print help message and exit
    if len(sys.argv)==1:
        psr.print_help(sys.stderr)
        sys.exit(1)


    args = psr.parse_args()
    args.func(args)

if __name__=="__main__":
    main()
