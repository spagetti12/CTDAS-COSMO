.DEFAULT_GOAL: clean

.PHONY: help clean

help:
	@echo " "
	@echo " Usage:"
	@echo "   make clean          # remove temporary (log) files"
	@echo " "

clean: 

# -------------------
	/bin/rm -f jb.[0-9]*.{jb,rc,log}
	/bin/rm -f *.pyc
	/bin/rm -f *.log
	/bin/rm -f *.[eo]*
	/bin/rm -f *.p[eo]*
	/bin/rm -f out.*
	/bin/rm -f err.*



