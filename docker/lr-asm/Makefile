TAG1 = kgarimella/lr-asm:0.01.09
TAG2 = kgarimella/lr-asm:latest

all: build push

build:
	docker build -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
