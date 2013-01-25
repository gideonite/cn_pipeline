LOGS_DIR = "logs"
NORMALIZED_DIR = "Level_2"
CBS_IN_DIR = "cbs_in"
CBS_OUT_DIR = "cbs_out"
GISTIC_IN_DIR = "gistic_in"
GISTIC_OUT_DIR = "gistic_in"


cbs:
	ln -s $(CBS_IN_DIR) cbs_in

all: $(NAME) tags

$(NAME): $(NAME).o
	gcc $(CFLAGS) $(NAME).o -o $(NAME)

$(NAME).o: $(NAME).c $(NAME).h
	gcc -c $(NAME).c

tags:
	ctags $(NAME).c $(NAME).h

clean:
	rm -f $(NAME) $(NAME).o
