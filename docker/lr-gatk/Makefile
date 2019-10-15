TAG1 = kgarimella/lr-gatk:0.01.00
TAG2 = kgarimella/lr-gatk:latest

all: build push

build:
	docker build -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
