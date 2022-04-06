VERSION = 0.5.1
TAG1 = quay.io/ymostovoy/modbam2bed:$(VERSION)
TAG2 = quay.io/ymostovoy/modbam2bed:latest

all: build push

build:
	docker build -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)