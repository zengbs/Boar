.PHONY: clean all

all: compile

compile:
	(cd src/integrator; $(MAKE) compile)
	(cd src/reconstruction; $(MAKE) compile)
	(cd src/rsolvers; $(MAKE) compile)
	(cd src; $(MAKE) compile)

clean:
	(cd src/integrator; $(MAKE) clean)
	(cd src/reconstruction; $(MAKE) clean)
	(cd src/rsolvers; $(MAKE) clean)
	(cd src; $(MAKE) clean)
