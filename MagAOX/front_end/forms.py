from flask_wtf import FlaskForm
from wtforms import StringField, FloatField, SubmitField
from wtforms.validators import InputRequired


class targetform(FlaskForm):
    targets = StringField('Targets', validators=[InputRequired()])
    starttime = StringField('Start Time', validators=[InputRequired()])
    endtime = StringField('End Time', validators=[InputRequired()])
    mode = StringField('Instrument Rotator Mode', validators=[InputRequired()])
    angle = FloatField('Instrument Rotator Angle', validators=[InputRequired()])
    guide1 = StringField('Guide Star 1', validators=[InputRequired()])
    guide2 = StringField('Guide Star 2', validators=[InputRequired()])
    pmepoch = FloatField('PM Epoch', validators=[InputRequired()])
    submit = SubmitField('Submit')                      
