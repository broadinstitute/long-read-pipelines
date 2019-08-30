TAG1 = kgarimella/lr-test:0.01.13
TAG2 = kgarimella/lr-test:latest

all: build push

build:
	docker build -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
