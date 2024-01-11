VERSION = 4.1.6
TAG1 = quay.io/ymostovoy/lr-repeatmasker:$(VERSION)
TAG2 = quay.io/ymostovoy/lr-repeatmasker:latest

all: build push

build:
	docker build -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1) 
	docker push $(TAG2)
