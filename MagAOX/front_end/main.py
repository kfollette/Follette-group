from flask import Flask, render_template, url_for
from forms import targetform

app = Flask(__name__)
app.config['SECRET_KEY'] = 'magaox_follettelab'

@app.route("/", methods=['GET', 'POST'])
def loadtargetform():
    return render_template('targetform.html', form=targetform())

if __name__ == '__main__':
    app.run(debug=True)
