TAG1 = kgarimella/lr-beryl:0.01.00
TAG2 = kgarimella/lr-beryl:latest

all: build push

build:
	docker build -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
