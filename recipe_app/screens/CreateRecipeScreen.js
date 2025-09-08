import React, { useState } from 'react';
import { View, TextInput, Button, Text, StyleSheet } from 'react-native';
import { auth, db } from '../services/firebaseConfig';
import { collection, addDoc } from 'firebase/firestore';

export default function CreateRecipeScreen({ navigation }) {
  const [title, setTitle] = useState('');
  const [ingredients, setIngredients] = useState('');
  const [steps, setSteps] = useState('');
  const [error, setError] = useState('');

  const submit = async () => {
    try {
      const ingList = ingredients.split(/,|\n/).map((i) => i.trim().toLowerCase()).filter(Boolean);
      await addDoc(collection(db, 'recipes'), {
        title,
        ingredients: ingList,
        steps,
        authorId: auth.currentUser.uid,
        likes: 0,
        likedBy: [],
        ratingTotal: 0,
        ratingCount: 0,
      });
      navigation.goBack();
    } catch (e) {
      setError(e.message);
    }
  };

  return (
    <View style={styles.container}>
      <Text style={styles.header}>New Recipe</Text>
      <TextInput style={styles.input} placeholder="Title" value={title} onChangeText={setTitle} />
      <TextInput style={styles.input} placeholder="Ingredients (comma separated)" value={ingredients} onChangeText={setIngredients} multiline />
      <TextInput style={styles.input} placeholder="Steps" value={steps} onChangeText={setSteps} multiline />
      {error ? <Text style={styles.error}>{error}</Text> : null}
      <Button title="Save" onPress={submit} color="#0066cc" />
    </View>
  );
}

const styles = StyleSheet.create({
  container: { flex: 1, padding: 20, backgroundColor: 'white' },
  header: { fontSize: 24, marginBottom: 20, color: '#0066cc' },
  input: { borderWidth: 1, borderColor: '#ddd', padding: 10, marginBottom: 10, borderRadius: 4 },
  error: { color: 'red', marginBottom: 10 }
});
