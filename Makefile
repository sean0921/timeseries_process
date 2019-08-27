#### Simple Makefile for build a pack

all:
	for i in $$(ls -d 0*/);do cd $$i;$(MAKE);cd -;done

mingw:
	for i in $$(ls -d 0*/);do cd $$i;$(MAKE) mingw;cd -;done

package:
	for i in $$(ls -d 0*/);do cd $$i;$(MAKE) package;cd -;done

mingw_package:
	for i in $$(ls -d 0*/);do cd $$i;$(MAKE) mingw_package;cd -;done

clean:
	for i in $$(ls -d 0*/);do cd $$i;$(MAKE) clean;cd -;done
