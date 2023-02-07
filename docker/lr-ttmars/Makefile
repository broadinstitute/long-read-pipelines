VERSION = 1.2
TAG1 = quay.io/ymostovoy/lr-ttmars:$(VERSION)
TAG2 = quay.io/ymostovoy/lr-ttmars:latest

all: build push

build:
	docker build -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)