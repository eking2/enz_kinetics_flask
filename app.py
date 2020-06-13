from flask import Flask, render_template
from flask_wtf import FlaskForm
from wtforms import SubmitField, IntegerField, DecimalField, StringField, TextAreaField
from wtforms.validators import DataRequired
from pathlib import Path
from kinetics import *

app = Flask(__name__)
app.config['SECRET_KEY'] = 'asdfk1j4kfjasf1kradf'
app.config['TMP'] = Path(Path.cwd(), 'tmp')
# do not cache
app.config['SEND_FILE_MAX_AGE_DEFAULT'] = 0

class EnzymeAssay(FlaskForm):

    total_rxn_vol = IntegerField('Total reaction volume (uL)', validators = [DataRequired()])
    enz_rxn_vol = IntegerField('Enzyme volume in reaction (uL)', validators = [DataRequired()])
    stock_conc = DecimalField('Enzyme stock concentration (mg/mL)', validators = [DataRequired()])
    enz_mol_wt = DecimalField('Enzyme molecular weight (g/mol)', validators = [DataRequired()])
    prot_dilution = IntegerField('Dilution factor', validators = [DataRequired()])
    extinct = DecimalField('Extinction coefficient (mM<sup>-1</sup> cm<sup>-1</sup>)', default = 6.22, validators = [DataRequired()])
    title = StringField('Run title')
    assay_data = TextAreaField('Assay data', render_kw={'rows' : 14, 'cols' : 70}, validators = [DataRequired()])

    submit = SubmitField('Submit')

@app.after_request
def add_header(response):

    response.headers['X-UA-Compatible'] = 'IE=Edge,chrome=1'
    response.headers['Cache-Control'] = 'public, max-age=0'
    return response


@app.route('/', methods = ['GET', 'POST'])
def home():
    form = EnzymeAssay()

    if form.validate():

        rxn_vol = form.total_rxn_vol.data
        enz_vol = form.enz_rxn_vol.data
        conc = form.stock_conc.data
        mol_wt = form.enz_mol_wt.data
        dil = form.prot_dilution.data
        ext = form.extinct.data
        title = form.title.data
        assay = form.assay_data.data

        # convert data
        assay_df = make_assay_df(assay)
        enz_dict = make_enz_dict(rxn_vol, enz_vol, conc, mol_wt, dil, ext)
        
        calc = kinetics_calc(assay_df, **enz_dict)
        calc.plot_mm(title)

        return render_template('index.html', result=True, form=form)

    return render_template('index.html', result = None, form=form)


if __name__ == '__main__':

    app.run(debug=True)
