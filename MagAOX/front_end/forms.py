from flask_wtf import FlaskForm
from wtforms import StringField, SubmitField
from wtforms.validators import InputRequired


class targetform(FlaskForm):
    targets = StringField('Targets',validators=[InputRequired()])
    #need to validate entry - query check?
    submit = SubmitField('Submit')                      
