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
    
    psr_bench.add_argument('-cf', '--config_file',
                              type=str,
                              required=True,
                              help="Relative path of config file."
                            )
    psr_bench.set_defaults(func=kpfcc.driver.bench)

    ## subcommand of astroq: plot -- run the plotting suite
    psr_plot = subpsr.add_parser('plot', parents=[psr_parent],
                                  description='Conduct benchmark tests',
                                  prefix_chars='-'
                                  )
    psr_plot.set_defaults(func='kpfcc.driver.plot')

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

    # If no arguments are provided, print help message and exit
    if len(sys.argv)==1:
        psr.print_help(sys.stderr)
        sys.exit(1)


    args = psr.parse_args()
    args.func(args)

if __name__=="__main__":
    main()
