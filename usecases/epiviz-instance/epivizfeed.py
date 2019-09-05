from epivizfeedcompute.server import setup_app, start_app
import os

if __name__ == "__main__":

    app = setup_app(os.getcwd() + "/epiviz.json")
    app.start_app()
    