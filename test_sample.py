from kpfcc import driver
import argparse

class TestClass:
    def test_one(self):
        x = "this"
        assert "h" in x

    def test_two(self):
        x = "hello"
        assert hasattr(x, "upper")
        
    def test_helloworld(self):
        
        args = argparse.Namespace(request_file='examples/hello-world.json',
                                  config_file='examples/config_hello_world.ini')
                                  
        x = driver.schedule(args)
        #assert isinstance(x, int)