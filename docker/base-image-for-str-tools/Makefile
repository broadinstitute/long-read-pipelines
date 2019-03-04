TAG1 = weisburd/base-image-for-str-tools:v1
TAG2 = weisburd/base-image-for-str-tools:latest

all: build push

build:
	docker build -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
