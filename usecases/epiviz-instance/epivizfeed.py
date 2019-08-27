from epivizfeedcompute.server import setup_app, start_app
import os

if __name__ == "__main__":

    app = setup_app(os.getcwd() + "/nsegil.json")
    app.start_app()
    