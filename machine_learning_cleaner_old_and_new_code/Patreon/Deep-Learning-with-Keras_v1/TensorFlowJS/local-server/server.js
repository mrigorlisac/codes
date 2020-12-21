let express = require("express"); // import express
let app = express(); //create express application, middleware funs

app.use(function(req, res, next) {
    console.log(`${new Date()} - ${req.method} request for ${req.url}`); // loging info about the request to the terminal where express server is running
    next(); // pass control to the next handler which ->
});

app.use(express.static("../static")); // -> responds with serving any static files placed here

app.listen(1100, function() { // specify what port express should listen, port 81. function() is executed
    console.log("Serving static on 1100");
});

// middleware functions are: function(req, res, next)
// and express.static(...)

// app.use only called once when the server is started running
// middleware functions inside them are called each time a request is sent