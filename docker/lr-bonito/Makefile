VERSION = 2.5.0 # This version number should follow the Megalodon version number
TAG1 = quay.io/ymostovoy/lr-megalodon:$(VERSION)
TAG2 = quay.io/ymostovoy/lr-megalodon:latest

all: build push

build:
	docker build -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)