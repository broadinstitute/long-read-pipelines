TAG1 = weisburd/gangstr:v1.4
TAG2 = weisburd/gangstr:latest

all: build push

build:
	docker build -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
