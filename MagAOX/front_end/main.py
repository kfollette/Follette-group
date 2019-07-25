from flask import Flask, render_template, url_for, flash, redirect
from forms import targetform
import catalog_entry

app = Flask(__name__)
app.config['SECRET_KEY'] = 'magaox_follettelab'

@app.route("/", methods=['GET', 'POST'])
def loadtargetform():
    form = targetform()
    if form.validate_on_submit():
        flash(f'Loading target information for {form.targets.data}!', 'success')
        return redirect(url_for('loadlandingpage'))
    return render_template('targetform.html', title='Target Form', form=form)

@app.route("/landingpage")
def loadlandingpage():
    return render_template('landingpage.html', title='Results')


if __name__ == '__main__':
    app.run(debug=True)
