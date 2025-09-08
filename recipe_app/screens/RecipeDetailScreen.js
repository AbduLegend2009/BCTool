import React, { useState } from 'react';
import { View, Text, Button, StyleSheet } from 'react-native';
import { doc, updateDoc, increment } from 'firebase/firestore';
import { db } from '../services/firebaseConfig';

export default function RecipeDetailScreen({ route }) {
  const { recipe } = route.params;
  const [rating, setRating] = useState(0);

  const rate = async () => {
    const ref = doc(db, 'recipes', recipe.id);
    await updateDoc(ref, { ratingTotal: increment(rating), ratingCount: increment(1) });
  };

  const avg = recipe.ratingCount ? (recipe.ratingTotal / recipe.ratingCount).toFixed(1) : '0';

  return (
    <View style={styles.container}>
      <Text style={styles.title}>{recipe.title}</Text>
      <Text style={styles.body}>Ingredients: {recipe.ingredients.join(', ')}</Text>
      <Text style={styles.body}>Steps: {recipe.steps}</Text>
      <Text style={styles.body}>Likes: {recipe.likes}</Text>
      <Text style={styles.body}>Rating: {avg}</Text>
      <View style={styles.ratingRow}>
        {[1,2,3,4,5].map((n) => (
          <Button key={n} title={`${n}`} onPress={() => setRating(n)} color={rating===n? '#0066cc':'#ccc'} />
        ))}
      </View>
      <Button title="Submit Rating" onPress={rate} color="#0066cc" />
    </View>
  );
}

const styles = StyleSheet.create({
  container: { flex: 1, padding: 20, backgroundColor: 'white' },
  title: { fontSize: 24, marginBottom: 10, color: '#0066cc' },
  body: { fontSize: 16, marginBottom: 10, color: '#0066cc' },
  ratingRow: { flexDirection: 'row', justifyContent: 'space-between', marginBottom: 10 }
});
