class RequestSet(object):
    """
    Request Set

    Specifies set of request to send to schedule. May be saved and restored from json
    """
    def __init__(self, meta, strategy, observable):
        self.meta = meta
        # enforce order
        strategy_cols = 'id t_visit n_intra_min n_intra_max tau_intra n_inter_max tau_inter'.split()
        strategy = strategy[strategy_cols]
        self.strategy = strategy
        self.observable = observable

    def __str__(self, n=10):
        s = "# Request Set # \n"
        s += "## Meta Data ## \n"
        s += '\n'*2
        s += self.meta.to_string()
        s += '\n'*2
        s += "## Strategy ##\n"
        s += self.strategy.head(n).to_string()
        s += '\n'*2
        s += "## Observable ##\n"
        s += self.observable.head(n).to_string()
        return s

    def to_json(self, fn):
        """
        Save request set to json
        """
        data = {}
        data['meta'] = self.meta.to_dict()
        data['strategy'] = self.strategy.to_dict()
        data['observable'] = self.observable.to_dict()

        with open(fn, "w") as f:
            json.dump(data,f,indent=4) # Pretty-printed for readability

def read_json(fn):
    """
    Read request set from json
    """
    with open(fn, "r") as f:
        data = json.load(f)

    meta = pd.Series(data['meta'])
    strategy = pd.DataFrame(data['strategy'])
    observable = pd.DataFrame(data['observable'])
    rs = RequestSet(meta, strategy, observable)
    return rs
