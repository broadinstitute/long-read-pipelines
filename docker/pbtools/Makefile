TAG1 = kgarimella/pbtools:v14
TAG2 = kgarimella/pbtools:latest

all: build push

build:
	docker build -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
