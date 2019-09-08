TAG1 = kgarimella/lr-align:0.01.17
TAG2 = kgarimella/lr-align:latest

all: build push

build:
	docker build -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
