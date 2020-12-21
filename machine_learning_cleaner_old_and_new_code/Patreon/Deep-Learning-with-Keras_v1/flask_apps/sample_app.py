from flask import Flask 

# creating instance of the flask class, __name__ is the name of the applicaton's
# module. Use it like this, if we are using the single module like here (in flask
# documentation).
app = Flask(__name__) 

# decorator that tells flask what url the user needs to browse to in order for the
# function underneath to be called. '/sample' -> endpoint, when user browse to this
# the function running will be called.
@app.route('/sample') 
def running():
    return 'FLASK IS RUNNING!!!!!'
