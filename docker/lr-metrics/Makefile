TAG1 = kgarimella/lr-metrics:0.01.05
TAG2 = kgarimella/lr-metrics:latest

all: build push

build:
	docker build -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
