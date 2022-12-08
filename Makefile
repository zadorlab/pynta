test: 
	nosetests --nocapture --nologcapture --verbose --with-coverage --cover-inclusive --cover-erase --cover-html --cover-html-dir=testing/coverage --cover-package=pynta pynta
