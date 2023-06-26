IDIR =  inc/
SDIR =  src/
ODIR =	obj/
DDIR =	dep/


csrc =  $(shell find $(SDIR) -name "*.c")
ccsrc = $(shell find $(SDIR) -name "*.cpp")
exsrc = $(wildcard GL/*.c)
OBJ = 	$(subst $(SDIR), $(ODIR), $(csrc:.c=.o)) \
		$(subst $(SDIR), $(ODIR), $(ccsrc:.cpp=.o)) \
		$(subst ../, $(ODIR),$(exsrc:.c=.o))
DEP =	$(subst $(ODIR), $(DDIR), $(OBJ:.o=.d))


#preprocessor flags
CPPFLAGS = -g -DDEBUG -DGLEW_STATIC -IGL -Isrc/
#linker flags (directories)
#LDFLAGS
#linker libraries
LDLIBS =-lm -pthread



all: macro.out

#common.out: $(CMNOBJ)
#	$(CXX) -o $@ $^ $(LDLIBS)

macro.out: $(OBJ)
	$(CXX) -o $@ $^ $(LDLIBS)



-include $(DEP)							#include all dep files in the makefile

.PHONY: clean
clean:
	rm -f $(OBJ)
	rm -f *.out
	rm -f $(DEP)


%/ : 
	@mkdir -p $@

.SECONDEXPANSION:
$(DDIR)%.d: $(SDIR)%.cpp | $$(dir $$@)
	@$(CPP) $(CPPFLAGS) $< -MM -MT $(subst $(DDIR), $(ODIR), $(@:.d=.o)) >$@


$(ODIR)%.o : $(SDIR)%.c | $$(dir $$@)
	$(CC) $(CPPFLAGS) -c $< -o $@

$(ODIR)%.o : ../%.c | $$(dir $$@)
	$(CC) $(CPPFLAGS) -c $< -o $@

$(ODIR)%.o : $(SDIR)%.cpp | $$(dir $$@)
	$(CXX) $(CPPFLAGS) -c $< -o $@


