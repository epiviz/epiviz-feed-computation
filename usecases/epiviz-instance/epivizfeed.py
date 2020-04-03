from epivizfeedcompute.server import setup_app, start_app
import os


if __name__ == "__main__":
    app = setup_app(server="http://54.157.53.251/api/", file=os.getcwd() + "/epiviz.json")
    # # app.start_app()
    start_app()
    