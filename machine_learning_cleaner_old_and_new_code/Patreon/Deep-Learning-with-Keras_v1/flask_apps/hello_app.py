from flask import request
from flask import jsonify
from flask import Flask
from flask import render_template

app = Flask(__name__)

# methods=['POST'] sending data along with the request
@app.route('/hello',methods=['POST'])
def hello():
    if request.method == 'POST':
        message = request.get_json(force=True)
        name = message['name']
        response = {'greeting': 'Hello, ' + name + '!'}
        return jsonify(response)
    else:
        return render_template('hello.html')
    
# # methods=['POST'] sending data along with the request
# @app.route('/hello',methods=['GET', 'POST'])
# def hello():
#     # request gives ability to access the data sent, gives the data in .json, force
#     # parsing .json even when unsure of data type
#     message = request.get_json(force=True)
#     name = message['name'] # values associated with the name key
#     response = {
#         'greeting': 'Hello, ' + name + '!'
#     }
#     return jsonify(response) # converts python dictionary into .json