from agents import search_refiner

def test_refiner_passthrough():
    queries = ["ion mobility spectroscopy protein"]
    search_refiner.settings.refine_search = False
    assert search_refiner.improve(queries) == queries
