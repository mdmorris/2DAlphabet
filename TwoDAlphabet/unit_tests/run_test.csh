python run_MLfit.py unit_tests/configs/input_lin100kv10k.json
python run_MLfit.py unit_tests/configs/input_lin100kv10kQuad.json
python run_Limit.py unit_tests/configs/input_lin100kv10k.json unitTest/lin100kv10k/
python run_Stats.py -d unitTest/lin100kv10k/ -t 100 --gof
python run_Stats.py -d unitTest/lin100kv10k/ -t 100 --diagnosticsWithToys
python run_Stats.py -d unitTest/lin100kv10k/ -t 100 --biasStudy -a unitTest/lin100kv10kQuad/
python run_Stats.py -d unitTest/lin100kv10k/ -t 100 --ftest -a unitTest/lin100kv10kQuad/