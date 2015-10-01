#
# Simple script to run all the slug problems in seciton 11
# and associated analysis script
#


echo "#############################################"
echo "#############################################"
echo "#############################################"
echo "Running all tests... This will take a few min"
echo "#############################################"
echo "#############################################"
echo "#############################################"

test/run_sfhsampling.sh
test/run_sampling.sh
test/run_imfchoice.sh
test/run_example_cluster.sh
test/run_example_galaxy.sh
test/run_constsampl.sh
test/run_cmfchoice.sh
test/run_clfraction.sh
test/run_cldisrupt.sh
test/run_spectra.sh

echo "#############################################"
echo "#############################################"
echo "#############################################"
echo "Done with tests... Running analysis"
echo "#############################################"
echo "#############################################"
echo "#############################################"

python test/plot_sfhsampling.py
python test/plot_sampling.py
python test/plot_imfchoice.py
python test/plot_example_cluster.py
python test/plot_example_galaxy.py
python test/plot_constsampl.py
python test/plot_cmfchoice.py
python test/plot_clfraction.py
python test/plot_cldisrupt.py
python test/plot_spectra.py

echo "#############################################"
echo "#############################################"
echo "#############################################"
echo "All done!!!!! "
echo "#############################################"
echo "#############################################"
echo "#############################################"
