import React, { useState, useEffect } from 'react';
import { View, Text, Button, FlatList, TextInput, StyleSheet, Image, TouchableOpacity } from 'react-native';
import { Camera } from 'expo-camera';
import { identifyIngredients } from '../services/GeminiService';
import { db, auth } from '../services/firebaseConfig';
import { collection, query, where, onSnapshot, doc, updateDoc, increment, arrayUnion } from 'firebase/firestore';

export default function HomeScreen({ navigation }) {
  const [cameraPermission, requestPermission] = Camera.useCameraPermissions();
  const [camera, setCamera] = useState(null);
  const [photo, setPhoto] = useState(null);
  const [ingredients, setIngredients] = useState([]);
  const [recipes, setRecipes] = useState([]);
  const [customIng, setCustomIng] = useState('');

  useEffect(() => {
    if (ingredients.length === 0) { setRecipes([]); return; }
    const q = query(collection(db, 'recipes'), where('ingredients', 'array-contains-any', ingredients));
    const unsub = onSnapshot(q, (snap) => {
      const r = snap.docs.map((d) => ({ id: d.id, ...d.data() }));
      setRecipes(r);
    });
    return () => unsub();
  }, [ingredients]);

  const takePicture = async () => {
    if (camera) {
      const photoData = await camera.takePictureAsync();
      setPhoto(photoData.uri);
      const list = await identifyIngredients(photoData.uri);
      setIngredients(list);
    }
  };

  const addIngredient = () => {
    if (customIng) {
      setIngredients([...ingredients, customIng.toLowerCase()]);
      setCustomIng('');
    }
  };

  const removeIngredient = (ing) => {
    setIngredients(ingredients.filter((i) => i !== ing));
  };

  const toggleLike = async (recipe) => {
    const ref = doc(db, 'recipes', recipe.id);
    await updateDoc(ref, { likes: increment(1), likedBy: arrayUnion(auth.currentUser.uid) });
  };

  const rateRecipe = async (recipe, rating) => {
    const ref = doc(db, 'recipes', recipe.id);
    await updateDoc(ref, { ratingTotal: increment(rating), ratingCount: increment(1) });
  };

  if (!cameraPermission) return <View />;
  if (!cameraPermission.granted) {
    return (
      <View style={styles.container}>
        <Text>We need your permission to use the camera</Text>
        <Button onPress={requestPermission} title="grant" />
      </View>
    );
  }

  return (
    <View style={styles.container}>
      {!photo ? (
        <Camera style={styles.camera} ref={setCamera}>
          <View style={styles.buttonContainer}>
            <TouchableOpacity style={styles.snap} onPress={takePicture} />
          </View>
        </Camera>
      ) : (
        <Image source={{ uri: photo }} style={styles.preview} />
      )}
      <View style={styles.ingContainer}>
        <FlatList
          data={ingredients}
          keyExtractor={(item) => item}
          horizontal
          renderItem={({ item }) => (
            <TouchableOpacity style={styles.ing} onPress={() => removeIngredient(item)}>
              <Text style={{ color: 'white' }}>{item}</Text>
            </TouchableOpacity>
          )}
        />
        <View style={styles.addRow}>
          <TextInput value={customIng} onChangeText={setCustomIng} placeholder="Add ingredient" style={styles.input} />
          <Button title="Add" onPress={addIngredient} color="#0066cc" />
        </View>
      </View>
      <FlatList
        data={recipes}
        keyExtractor={(item) => item.id}
        renderItem={({ item }) => (
          <TouchableOpacity style={styles.recipe} onPress={() => navigation.navigate('RecipeDetail', { recipe: item })}>
            <Text style={styles.recipeTitle}>{item.title}</Text>
            <Button title="Like" onPress={() => toggleLike(item)} color="#0066cc" />
          </TouchableOpacity>
        )}
      />
      <Button title="Profile" onPress={() => navigation.navigate('Profile')} color="#0066cc" />
      <Button title="Create Recipe" onPress={() => navigation.navigate('CreateRecipe')} color="#0066cc" />
    </View>
  );
}

const styles = StyleSheet.create({
  container: { flex: 1, backgroundColor: 'white' },
  camera: { flex: 1 },
  buttonContainer: { flex: 1, backgroundColor: 'transparent', justifyContent: 'flex-end' },
  snap: { alignSelf: 'center', marginBottom: 20, width: 70, height: 70, borderRadius: 35, backgroundColor: '#0066cc' },
  preview: { flex: 1 },
  ingContainer: { padding: 10 },
  ing: { backgroundColor: '#0066cc', padding: 6, marginRight: 6, borderRadius: 4 },
  addRow: { flexDirection: 'row', alignItems: 'center', marginTop: 10 },
  input: { flex: 1, borderWidth: 1, borderColor: '#ddd', marginRight: 10, padding: 8 },
  recipe: { padding: 10, borderBottomWidth: 1, borderBottomColor: '#eee' },
  recipeTitle: { fontSize: 18, color: '#0066cc' }
});
