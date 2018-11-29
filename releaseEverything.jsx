 #target illustrator

// releaseEverything.jsx
//
// Copyright (c) 2017 Janne Ojala
//
// Licence: https://opensource.org/licenses/MIT

(function(){


//var objs = app.activeDocument.pageItems;
var objs = app.activeDocument.selection;
// or if you want only selection use app.activeDocument.selection
traverseSceneObjects(objs);


function traverseSceneObjects(pageItems){

    for (var iter=0 ; iter<pageItems.length; iter++ ){
        var item = pageItems[iter];
        var typename = item.typename;

        // apply action or get the subitems of object
        if (typename === "PathItem"){
            item.clipping = false;

        } else if (typename === "GroupItem") {
            traverseSceneObjects( item.pageItems );
            release( item, "pageItems" );

        } else if (typename === "CompoundPathItem" ) {
            traverseSceneObjects( item.pathItems );
            release( item, "pathItems" );
        }

    }

}


function release(obj, action) {

    for (var i=obj[action].length-1 ; i>=0; i--){
        obj[action][i].move( obj, ElementPlacement.PLACEAFTER );

    }

}

})();