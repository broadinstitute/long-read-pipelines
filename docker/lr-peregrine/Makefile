TAG1 = kgarimella/lr-peregrine:0.01.00
TAG2 = kgarimella/lr-peregrine:latest

all: build push

build:
	docker build -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
