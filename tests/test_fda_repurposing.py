"""Scientific data-contract regression tests; synthetic records only."""
import importlib.util
import contextlib
import csv
import io
import json
from pathlib import Path
import sqlite3
import tempfile
import unittest
import zipfile

try:
    import pandas as pd
except ImportError:
    pd = None

if pd is not None:
    spec = importlib.util.spec_from_file_location('fda', Path(__file__).resolve().parents[1]/'pipeline/fda_repurposing.py')
    fda = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(fda)


@unittest.skipIf(pd is None, 'Run with the chemical environment (pandas/pyarrow)')
class TestFDA(unittest.TestCase):
    def tables(self):
        ids = [('001','1','ALPHA'), ('002','3','ALPHA SALT'),
               ('003','4','BETA'), ('004','1','ALPHA; BETA'),
               ('005','1','GAMMA'), ('006','1','BETA')]
        return {
            'Products': [{'ApplNo':a,'ProductNo':'001','ActiveIngredient':n,
                          'DrugName':n,'Form':'TABLET;ORAL','Strength':'1MG'} for a,s,n in ids],
            'Applications': [{'ApplNo':a,'ApplType':'NDA'} for a,s,n in ids],
            'Submissions': [{'ApplNo':a,'SubmissionType':'ORIG','SubmissionStatus':'AP',
                             'SubmissionStatusDate':'2000-01-01'} for a,s,n in ids if a!='005'],
            'MarketingStatus': [{'ApplNo':a,'ProductNo':'001','MarketingStatusID':s} for a,s,n in ids],
            'MarketingStatus_Lookup': [{'MarketingStatusID':s,'MarketingStatusDescription':label}
                    for s,label in [('1','Prescription'),('3','Discontinued'),('4','None (Tentative Approval)')]]}

    def products(self):
        return fda.fda_products(self.tables())

    def database(self):
        c = sqlite3.connect(':memory:')
        c.executescript('''
          CREATE TABLE molecule_dictionary (molregno INTEGER PRIMARY KEY, chembl_id TEXT, pref_name TEXT);
          INSERT INTO molecule_dictionary VALUES (1,'CHEMBL1','ALPHA'),(2,'CHEMBL2','ALPHA SALT'),(3,'CHEMBL3','BETA');
          CREATE TABLE molecule_synonyms (molregno INTEGER, synonyms TEXT);
          INSERT INTO molecule_synonyms VALUES (1,'AMBIGUOUS'),(3,'AMBIGUOUS');
          CREATE TABLE molecule_hierarchy (molregno INTEGER, parent_molregno INTEGER);
          INSERT INTO molecule_hierarchy VALUES (1,1),(2,1),(3,3);
          CREATE TABLE compound_structures (molregno INTEGER, canonical_smiles TEXT, standard_inchi_key TEXT);
          INSERT INTO compound_structures VALUES (1,'CC','KEY1'),(3,'CCC','KEY3');
          CREATE TABLE target_dictionary (tid INTEGER, chembl_id TEXT, pref_name TEXT, organism TEXT, target_type TEXT);
          INSERT INTO target_dictionary VALUES (1,'T1','target','Homo sapiens','SINGLE PROTEIN');
          CREATE TABLE target_components (tid INTEGER, component_id INTEGER);
          INSERT INTO target_components VALUES (1,1);
          CREATE TABLE component_sequences (component_id INTEGER, accession TEXT);
          INSERT INTO component_sequences VALUES (1,'P1');
          CREATE TABLE component_synonyms (component_id INTEGER, component_synonym TEXT, syn_type TEXT);
          INSERT INTO component_synonyms VALUES (1,'GENE1','GENE_SYMBOL');
          CREATE TABLE assays (assay_id INTEGER, tid INTEGER, chembl_id TEXT, assay_type TEXT,
            description TEXT, assay_organism TEXT, assay_cell_type TEXT, assay_tissue TEXT, confidence_score INTEGER);
          INSERT INTO assays VALUES (1,1,'A1','B','binding','Homo sapiens',NULL,NULL,9);
          CREATE TABLE activities (activity_id INTEGER, assay_id INTEGER, molregno INTEGER, doc_id INTEGER,
            standard_type TEXT, standard_relation TEXT, standard_value REAL, standard_units TEXT,
            pchembl_value REAL, activity_comment TEXT, data_validity_comment TEXT, potential_duplicate INTEGER,
            action_type TEXT);
          INSERT INTO activities VALUES (1,1,1,1,'IC50','=',10,'nM',8,NULL,NULL,0,NULL),
             (2,1,2,1,'IC50','>',10000,'nM',NULL,NULL,NULL,0,NULL);
          CREATE TABLE docs (doc_id INTEGER, chembl_id TEXT, pubmed_id TEXT, doi TEXT, year INTEGER, title TEXT);
          INSERT INTO docs VALUES (1,'D1','123','doi',2001,'source');
        ''')
        self.addCleanup(c.close)
        return c

    def test_tentative_and_missing_original_approval_excluded(self):
        p = self.products()
        self.assertFalse(p[2]['approval_evidence_confirmed'])
        self.assertFalse(p[4]['eligible_current_product'])
        self.assertTrue(p[0]['eligible_current_product'])

    def test_discontinued_not_equated_with_safety_withdrawal(self):
        p = self.products()[1]
        self.assertTrue(p['approval_evidence_confirmed'])
        self.assertFalse(p['eligible_current_product'])
        self.assertEqual(p['safety_withdrawal_status'], 'not_assessed')

    def test_exact_mapping_quarantines_ambiguity_and_collapses_salt(self):
        r, issues = fda.resolve_names(self.database(), {'ALPHA','ALPHA SALT','AMBIGUOUS','ALPH'})
        self.assertEqual(r['ALPHA']['ingredient_id'], r['ALPHA SALT']['ingredient_id'])
        self.assertEqual({i['issue'] for i in issues}, {'ambiguous_parent','unmatched'})

    def test_combination_preserved_even_with_unresolved_component(self):
        r,_ = fda.resolve_names(self.database(), {'ALPHA','ALPHA SALT'})
        products = self.products()
        ingredients = fda.build_catalog(products, r)
        self.assertEqual(products[3]['mapped_ingredient_ids'], ['CHEMBL1',None])
        self.assertTrue(products[3]['is_combination'])
        self.assertEqual(ingredients[0]['eligible_single_ingredient_product_ids'], ['001:001'])

    def test_measured_conflicts_censoring_and_unknown_direction_retained(self):
        c = self.database()
        r,_ = fda.resolve_names(c, {'ALPHA','ALPHA SALT','BETA'})
        ingredients = fda.build_catalog(self.products(), r)
        frame, coverage = fda.measured_interactions(c, ingredients, ['GENE1','MISSING'])
        self.assertEqual(len(frame), 2)
        self.assertEqual(set(frame.standard_relation), {'=','>'})
        self.assertEqual(set(frame.tested_molecule_chembl_id), {'CHEMBL1','CHEMBL2'})
        self.assertEqual(set(frame.direction_status), {'unknown'})
        self.assertEqual(set(coverage.gene), {'GENE1'})

    def test_sibling_join_by_parent_and_no_invented_endpoint(self):
        c = self.database()
        r,_ = fda.resolve_names(c, {'ALPHA'})
        ingredients = fda.build_catalog(self.products(), r)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'siblings.parquet'
            pd.DataFrame([{'chembl_id':'CHEMBL2','disease':'Syndrome A','intent':'instability_directed'},
                          {'chembl_id':'CHEMBL3','disease':'Syndrome A','intent':'consequence_directed'}]).to_parquet(path)
            result = fda.sibling_candidates(c, ingredients, path, ['Syndrome A'])
        self.assertEqual(len(result),1)
        self.assertEqual(result.iloc[0].ingredient_id,'CHEMBL1')
        self.assertEqual(result.iloc[0].measured_endpoint,'not_retrieved')

    def test_nan_direction_is_unknown_in_mixed_realistic_column(self):
        c = self.database()
        c.execute("UPDATE activities SET action_type='INHIBITOR' WHERE activity_id=1")
        r,_ = fda.resolve_names(c, {'ALPHA'})
        frame,_ = fda.measured_interactions(c, fda.build_catalog(self.products(),r), ['GENE1'])
        self.assertEqual(frame.set_index('activity_id').loc[2,'direction_status'],'unknown')

    def test_retrieval_cannot_nominate(self):
        r = fda.candidate_record({'canonical_name':'example'},['related_disease'],[])
        self.assertEqual(r['decision'],'insufficient_evidence')
        self.assertEqual(r['assessments']['functional_rescue'],'unestablished')

    def test_full_run_provenance_historical_results_and_no_overwrite(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            with zipfile.ZipFile(root/'fda.zip','w') as z:
                for table, rows in self.tables().items():
                    buf = io.StringIO()
                    writer = csv.DictWriter(buf, fieldnames=list(rows[0]), delimiter='\t')
                    writer.writeheader()
                    writer.writerows(rows)
                    z.writestr(table+'.txt',buf.getvalue())
            con = sqlite3.connect(root/'chembl.db')
            self.database().backup(con)
            con.close()
            pd.DataFrame([{'chembl_id':'CHEMBL2','disease':'Syndrome A','intent':'instability_directed'}]).to_parquet(root/'sibling.parquet')
            config = {'fda_zip':str(root/'fda.zip'), 'chembl_db':str(root/'chembl.db'),
                      'sibling_table':str(root/'sibling.parquet'), 'genes':['GENE1'],
                      'sibling_diseases':['Syndrome A']}
            out = root/'out'
            with contextlib.redirect_stdout(io.StringIO()):
                fda.run(config,out)
            summary = json.loads((out/'summary.json').read_text())
            self.assertEqual(summary['candidate_ingredients'],1)
            self.assertEqual(summary['measured_records'],2)
            manifest = json.loads((out/'manifest.json').read_text())
            self.assertEqual(manifest['output_hashes']['candidates.json'], fda.sha256(out/'candidates.json'))
            self.assertEqual(manifest['implementation_sha256'],fda.sha256(fda.__file__))
            historical = json.loads((out/'baseline_candidates.json').read_text())
            record = json.loads((out/'baseline.json').read_text())
            # Recorded run artifacts are retained privately; either state is provenance, never a nomination.
            if record['historical_baseline'] == 'available':
                self.assertTrue(historical)
                for r in historical:
                    self.assertEqual(r['original_result']['gate']['status'],'review_required')
            else:
                self.assertEqual(historical,[])
                self.assertIn('not distributed',record['reason'])
            with self.assertRaisesRegex(ValueError,'fresh'):
                fda.run(config,out)

    def test_historical_baseline_preserves_recorded_gates_when_available(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source, out = root/'source', root/'out'
            source.mkdir(); out.mkdir()
            (source/'04_candidate_01.json').write_text(json.dumps(
                {'candidate':{'name':'Molecule A','identity':{'chembl_id':'CHEMBL1'},'evidence':[]},
                 'gate':{'status':'review_required'}}))
            fda.baseline(out, source)
            historical = json.loads((out/'baseline_candidates.json').read_text())
            self.assertEqual(len(historical),1)
            self.assertEqual(historical[0]['original_result']['gate']['status'],'review_required')
            self.assertEqual(historical[0]['decision_reason'],
                             'Historical result preserved; no FDA eligibility inferred')
            self.assertEqual(json.loads((out/'baseline.json').read_text())['historical_baseline'],'available')


if __name__ == '__main__':
    unittest.main()
