VERSION = 0.1.9
TAG1 = quay.io/ymostovoy/lr-xtea:$(VERSION)
TAG2 = quay.io/ymostovoy/lr-xtea:latest

all: build push

build:
	docker build -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
