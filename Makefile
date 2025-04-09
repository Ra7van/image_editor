CC=gcc
CFLAGS=-Wall -Wextra -std=c99
TARGETS=image_editor

build: $(TARGETS)

image_editor: main.c functii.c
	$(CC) $(CFLAGS) -g -o image_editor main.c functii.c -lm

pack:
	zip -FSr 3XYCA_FirstnameLastname_Tema3.zip README Makefile *.c *.h

clean:
	rm -f $(TARGETS)

.PHONY: pack clean
