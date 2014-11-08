#
# Simple script to run the slug problem: example, then
# use cloudy to compute the corresponding nebular emission. 
# See section 10 of the SLUG manual. 
#

python $SLUG_DIR/bin/slug.py $SLUG_DIR/param/example_cloudy.param 
python $SLUG_DIR/cloudy_slug/cloudy_slug.py $SLUG_DIR/output/SLUG_CLOUDY -v
