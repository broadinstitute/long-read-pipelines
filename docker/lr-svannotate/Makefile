VERSION = 4.5.0.0
TAG1 = quay.io/ymostovoy/lr-svannotate:$(VERSION)
TAG2 = quay.io/ymostovoy/lr-svannotate:latest

all: build push

build:
	docker build -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
