import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'pipeline'))
import run
from prepare_public import Paragraphs
sys.path.insert(0, str(ROOT/'mva_llm/src'))
from response_parser import parse_calls, user_answer


class PortablePipeline(unittest.TestCase):
    def config(self):
        return json.loads((ROOT/'config/example.json').read_text())

    def test_all_stages_and_existing_implementations(self):
        plan = run.build_plan(self.config(), run.STAGES)
        self.assertEqual(list(dict.fromkeys(x['stage'] for x in plan)), list(run.STAGES))
        for item in plan:
            for argument in item['argv']:
                if argument.endswith(('.py','.sh')):
                    self.assertTrue(Path(argument).is_file(), argument)
        commands = [Path(x['argv'][1]).name for x in plan if len(x['argv']) > 1]
        self.assertIn('train_mlp.py', commands)
        self.assertIn('train_targets.py', commands)
        self.assertLess(commands.index('check_resources.py'), commands.index('workflow.py'))

    def test_relocated_checkout_with_spaces(self):
        previous = run.ROOT
        with tempfile.TemporaryDirectory(prefix='checkout with spaces ') as d:
            run.ROOT = Path(d)
            try:
                c = self.config()
                c['patient_vcf'] = 'inputs/a file.vcf.gz'
                plan = run.build_plan(c, ['patient','reason'])
                prep = next(x for x in plan if x['argv'][1].endswith('prepare_patient.py'))
                self.assertIn(str(Path(d)/'inputs/a file.vcf.gz'), prep['argv'])
                generated = run.workflow_config(c)
                self.assertTrue(generated['variants'].startswith(d))
                self.assertNotIn('NMN', json.dumps(generated))
            finally: run.ROOT = previous

    def test_reference_drift_requires_explicit_option(self):
        c = self.config()
        plan = run.build_plan(c, ['reason'])
        self.assertNotIn('--allow-reference-drift', plan[0]['argv'])
        c['allow_reference_drift'] = True
        self.assertIn('--allow-reference-drift', run.build_plan(c, ['reason'])[0]['argv'])

    def test_model_channels_do_not_leak_reasoning(self):
        self.assertEqual(user_answer('to=self private to=user answer to=self more to=user end'), 'answer\n\nend')
        self.assertIsNone(user_answer('to=self unfinished'))
        call = '<atem:invoke name="paper_text"><atem:parameter name="pmid">123</atem:parameter></atem:invoke>'
        self.assertEqual(parse_calls(call), [{'name':'paper_text','arguments':{'pmid':123}}])

    def test_paper_text_preserves_scientific_notation(self):
        parser = Paragraphs()
        parser.feed('<p>NAD<sup>+</sup> and BubR1<sup>H/H</sup>\n  cells.</p><p></p>')
        self.assertEqual(parser.paragraphs, ['NAD+ and BubR1H/H cells.', ''])


if __name__ == '__main__': unittest.main()
