test-all:
	nosetests --nocapture --nologcapture --verbose --with-coverage --cover-inclusive --cover-erase --cover-html --cover-html-dir=testing/coverage --cover-package=pynta pynta

test-unittests:
	nosetests --nocapture --nologcapture -A 'not functional' --verbose --with-coverage --cover-inclusive --cover-erase --cover-html --cover-html-dir=testing/coverage --cover-package=pynta pynta

test-functional:
	nosetests --nocapture --nologcapture -A 'functional' --verbose --with-coverage --cover-inclusive --cover-erase --cover-html --cover-html-dir=testing/coverage --cover-package=pynta pynta
