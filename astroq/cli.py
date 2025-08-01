"""
Command Line Interface
"""

# Standard library imports
import argparse
import os
import sys
from configparser import ConfigParser

# Third-party imports
import numpy as np
import pandas as pd

# Local imports
import astroq

# Create the parser at module level
parser = argparse.ArgumentParser(
    description='AstroQ: Optimized observation scheduling',
    prog='astroq'
)

def main():
    # This parser defines the top-level command "astroq"
    psr = parser  # Use the global parser

    psr.add_argument('-V', '--version',
                     action='version',
                     version="%(prog)s {}".format(astroq.__version__),
                     help="Print version number and exit."
                     )
    # Parsers spawned from this subparser generator will define subcommands of the astroq top-level command
    # Parsers spawned from subpsr (so it defines a subcommand) has psr_parent as its parent (so it accepts any arguments of psr_parent in addition to its own)
    subpsr = psr.add_subparsers(title='subcommands', dest='subcommand', required=True)
    ## Parent parser to define arguments common to all subcommands (so you don't have to redefine them)
    psr_parent = argparse.ArgumentParser(add_help=False)

    ## subcommand of astroq: bench -- conduct benchmark tests
    psr_bench = subpsr.add_parser('bench', parents=[psr_parent],
                                  description='Conduct benchmark tests',
                                  prefix_chars='-'
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
    psr_bench.add_argument('-t', '--thin',
                              type=int,
                              default=1,
                              help="Downsample the request frame by this factor for faster testing (default: 1, no thinning)."
                            )
    psr_bench.set_defaults(func=astroq.driver.bench)

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
    psr_plot.set_defaults(func=astroq.driver.plot_pkl)


    ## subcommand of astroq: kpfcc -- Do KPFCC stuff
    psr_kpfcc = subpsr.add_parser('kpfcc', parents=[psr_parent],
                                  description='Do KPFCC stuff',
                                  prefix_chars='-'
                                  )
    kpfcc_subpsr = psr_kpfcc.add_subparsers(title='kpfcc subcommands', dest='kpfcc_subcommand')
    psr_kpfcc.set_defaults(func=astroq.driver.kpfcc)

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
    psr_kpfcc_prep.add_argument('-as', '--allo_source',
                            type=str,
                            default="db",
                            help="Absolute path of observatory-provided allocation file. Use 'db' to pull from the database."
                              )
    psr_kpfcc_prep.add_argument('-ps', '--past_source',
                            type=str,
                            default="db",
                            help="Absolute path of a past history file. Use 'db' to pull from the database."
                              )
    psr_kpfcc_prep.set_defaults(func=astroq.driver.kpfcc_prep)

 
    ## subcommand of kpfcc: webapp -- launch web app to view interactive plots
    psr_kpfcc_webapp = kpfcc_subpsr.add_parser('webapp', #parents=[psr_parent],
                                               description="Launch web app to view interactive plots",
                                               prefix_chars="-"
                                               )
    psr_kpfcc_webapp.add_argument('-cf', '--config_file',
                              type=str,
                              required=True,
                              help="Path to config file."
                                )
    psr_kpfcc_webapp.set_defaults(func=astroq.driver.kpfcc_webapp)

    ## subcommand of kpfcc: plan-semester -- plan a semester's worth of observations
    psr_kpfcc_plan_semester = kpfcc_subpsr.add_parser('plan-semester', #parents=[psr_parent],
                                                       description="Plan a semester's worth of observations using optimization (builds request set on-the-fly)",
                                                       prefix_chars="-"
                                                       )
    psr_kpfcc_plan_semester.add_argument('-cf', '--config_file',
                                         type=str,
                                         required=True,
                                         help="Relative path of config file."
                                         )
    psr_kpfcc_plan_semester.set_defaults(func=astroq.driver.kpfcc_plan_semester)

    ## subcommand of kpfcc: plan-night -- run the night planner
    psr_kpfcc_plan_night = kpfcc_subpsr.add_parser('plan-night', #parents=[psr_parent],
                                                    description="Run the night planner (Target & Time Planner)",
                                                    prefix_chars="-"
                                                    )
    psr_kpfcc_plan_night.add_argument('-cf', '--config_file',
                                      type=str,
                                      required=True,
                                      help="Relative path of config file."
                                      )
    psr_kpfcc_plan_night.set_defaults(func=astroq.driver.ttp)

    ## subcommand of astroq: compare -- compare request set and schedule file
    psr_compare = subpsr.add_parser('compare', parents=[psr_parent],
                                    description='Compare request set and schedule for consistency',
                                    prefix_chars='-'
                                    )
    psr_compare.add_argument('-cf', '--config_file',
                              type=str,
                              required=True,
                              help="Relative path of config file."
                              )
    psr_compare.add_argument('-sf', '--schedule_file',
                              type=str,
                              required=True,
                              help="Relative path of schedule file."
                              )
    psr_compare.set_defaults(func=astroq.driver.requests_vs_schedule)

    ## subcommand of astroq: simsemester -- simulate a semester with a given weather loss pattern.
    # psr_simsemester = subpsr.add_parser('simsemester', parents=[psr_parent],
    #                                 description='Compare request set and schedule for consistency',
    #                                 prefix_chars='-'
    #                                 )
    # psr_simsemester.add_argument('-cf', '--config_file',
    #                           type=str,
    #                           required=True,
    #                           help="Relative path of config file."
    #                           )
    # psr_simsemester.set_defaults(func=astroq.driver.make_simulated_history)

    # If no arguments are provided, print help message and exit
    if len(sys.argv)==1:
        psr.print_help(sys.stderr)
        sys.exit(1)

    args = psr.parse_args()
    args.func(args)

if __name__=="__main__":
    main()
