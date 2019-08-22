#### Simple Makefile for build a pack

all:
	for i in $$(ls -d 0*/);do cd $$i;$(MAKE);cd -;done

mingw:
	for i in $$(ls -d 0*/);do cd $$i;$(MAKE) mingw;cd -;done

install:
	for i in $$(ls -d 0*/);do cd $$i;$(MAKE) install;cd -;done

mingw_install:
	for i in $$(ls -d 0*/);do cd $$i;$(MAKE) mingw_install;cd -;done

clean:
	for i in $$(ls -d 0*/);do cd $$i;$(MAKE) clean;cd -;done
