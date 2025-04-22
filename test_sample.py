from kpfcc import driver
import argparse

class TestClass:

    def test_bench(self):
        # override the tests of "do these files exist" and just run it all
        # but run a shortened version, so it doesn't take long, see shortcut=5
        bn.build_toy_model_from_paper(savepath="examples/bench/inputs/",shortcut=5)
        dr.kpfcc_prep(argparse.Namespace(config_file='examples/bench/config_benchmark.ini'))
        dr.kpfcc_build(argparse.Namespace(config_file='examples/bench/config_benchmark.ini'))
        bn.do_benchmark_files_exist('examples/bench/config_benchmark.ini', shortcut=0)
        request_set = rq.read_json("examples/bench/outputs/2024-08-01/request_set.json")
        request_set = bn.firstN_Requests(2, request_set)
        request_set = bn.set_nSlots_singles(2, request_set, start_row=1)
        request_set.to_json("examples/bench/outputs/2024-08-01/request_set_short.json")
        dr.schedule(argparse.Namespace(request_file="examples/bench/outputs/2024-08-01/request_set_short.json", config_file='examples/bench/config_benchmark.ini'))             

    def test_prep(self):
        # test the creation of all upstream files, including the allocation map 
        dr.kpfcc_prep(argparse.Namespace(config_file='examples/hello_world/config_hello_world.ini'))
