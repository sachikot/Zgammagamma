CPPFLAGS=`root-config --cflags` -I/usr/local/root/include
LDFLAGS=`root-config --libs` -lMinuit
skimmer_eeg: skimmer_eeg.cc
	g++ $(CPPFLAGS) $(LDFLAGS) -o $@ $^
