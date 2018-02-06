// google.maps.event.addDomListener(window, 'load', () => alert('hi'));
// let decoded = google.maps.geometry.encoding.decodePath('kh~uElvolR|C`@`I|@nEh@lEh@rC\v@HjCZ`AD~ANhBR|BVxFr…\`@DdALfCXnEh@f@FdD`@|ARjBRlEf@tBTxARjEf@tC\LBj@H');
// console.log(decoded);

var getClosestPoint = function (location, path, toleranceInMeters) {
    let closestPoint = null;
    let minDistance = 99999999999;
    for (var currentPoint of path) {
        let currentDistance = google.maps.geometry.spherical.computeDistanceBetween(location, currentPoint);
        if (currentDistance <= toleranceInMeters && currentDistance < minDistance) {
            minDistance = currentDistance;
            closestPoint = currentPoint;
            // return true;
        }
    }

    return closestPoint;
    // return false;
};

window.addEventListener('load', initializeMap);

function initializeMap() {
    let options = {
        center: new google.maps.LatLng(42.65494, 23.336084),
        zoom: 16
    };

    let domContainer = document.getElementById('map');
    let map = new google.maps.Map(domContainer, options);

    let directionsService = new google.maps.DirectionsService;
    let directionsDisplay = new google.maps.DirectionsRenderer;
    directionsDisplay.setMap(map);

    let passangerDirectionsService = new google.maps.DirectionsService;
    let passangerDirectionsDisplay = new google.maps.DirectionsRenderer;
    directionsDisplay.setMap(map);

    let onChangeHandler = () => {
        let startEl = document.getElementById('start');
        let endEl = document.getElementById('end');
        if (!startEl.value || !endEl.value) {
            return;
        }

        calculateAndDisplayRoute(startEl.value, endEl.value, directionsService, directionsDisplay, 'DRIVING', map);
    };

    document.getElementById('btn-display').addEventListener('click', onChangeHandler);

    let marker = new google.maps.Marker({ map });
    // let poly = new google.maps.Polyline({
    //     map: map,
    //     path: [
    //         new google.maps.LatLng(42.65666478, 23.33636212),
    //         new google.maps.LatLng(42.65378476, 23.33606171),
    //     ],
    //     strokeColor: '#FF0000',
    //     strokeOpacity: 1.0,
    //     strokeWeight: 2
    // });

    map.addListener('click', (e) => {
        let newPos = {
            lat: e.latLng.lat(),
            lng: e.latLng.lng()
        };

        marker.latLng
        marker.setPosition(newPos);

        let routePath = route.overview_path;
        let toleranceMeter = 1000;
        let closestPoint = getClosestPoint(e.latLng, route.overview_path, toleranceMeter);
        if (closestPoint == null) {
            alert(`No results found.\r\nTry range wider than ${toleranceMeter}`);
            return;
        }

        console.log(route.overview_path);
        let walkingPolyLine = new google.maps.Polyline({
            map: map,
            path: markers.map(m => m.position),
            strokeColor: '#FF0000',
            strokeOpacity: 1.0,
            strokeWeight: 2
        });

        calculateAndDisplayRoute(
            newPos,
            closestPoint,
            passangerDirectionsService,
            passangerDirectionsDisplay,
            'WALKING',
            map
        );

        // let tolerance = 0.01;
        // if (google.maps.geometry.poly.isLocationOnEdge(marker.position, pathPolyLine, tolerance)) {
        //     console.log(`You are within range with ${tolerance} tolerance.`)
        // }
    });
}

var route;
function calculateAndDisplayRoute(start, end, directionsService, directionsDisplay, travelMode, map) {
    directionsService.route({
        origin: start,
        destination: end,
        travelMode: travelMode,
    }, (response, status) => {
        if (status === 'OK') {
            // if (travelMode === 'DRIVING') {
                route = response.routes[0];
            // }

            let count = 0;
            response.routes[0].legs[0].steps.forEach(s => count += s.path.length);
            console.log('Points count: ', count);
            console.log('Overview Points count: ', response.routes[0].overview_path.length);
            console.log(response);
            setMarkers(response.routes[0].overview_path, map);
            // response.routes[0].legs[0].steps.forEach(step => setMarkers(step.path, map))
            directionsDisplay.setDirections(response);
        } else {
            window.alert('Directions request failed due to ' + status);
        }
    });
}

var markers;
function setMarkers(path, map) {
    markers = [];
    path.forEach(e => {
        let marker = new google.maps.Marker({
            position: {
                lat: e.lat(),
                lng: e.lng()
            },
            map: map
        });

        markers.push(marker);
    });

    console.log(markers);
}