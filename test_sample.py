class TestClass:
    def test_one(self):
        x = "this"
        assert "h" in x

    def test_two(self):
        x = "hello"
        assert hasattr(x, "upper")
        
    def test_helloworld(self):
        x = some_output_of_hello_world
        assert hasattr(x, 'some_attribute')